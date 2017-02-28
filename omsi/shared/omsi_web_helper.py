"""
Module with helper functions for interactions with the OpenMSI web infrastructure,
e.g. update job status, explicitly add a file to the OpenMSI database,
update file permissions so that Apache can access it etc.
"""
import os
import warnings
import sys
import time
# Imports for user input with timeout
from select import select
import platform
try:
    if platform.system() == "Windows":
        import msvcrt
except ImportError:
    pass
from omsi.shared.log import log_helper


class WebHelper(object):
    """
    Class providing a collection of functions for web-related file conversion
    tasks, e.g, : i) adding files to the web database, ii) notifying users via email,
    iii) setting file permissions for web-access.
    """

    default_db_server_url = "https://openmsi.nersc.gov/"
    allowed_nersc_locations = ["/project/projectdirs/openmsi/omsi_data_private",
                               "/global/project/projectdirs/openmsi/omsi_data_private",
                               "/data/openmsi/omsi_data"]
    super_users = ['bpb', 'oruebel']  # Simple safe-guard to prevent non-admins to add files from arbitrary locations

    def __init__(self):
        pass

    @staticmethod
    def upload_file_to_nersc(filepath,
                             username=None,
                             target_dir=None,
                             session=None,
                             persist_session=False,
                             machine='cori',
                             register=False):
        """
        Upload a file to NERSC via NEWT

        :param filepath: The path of the file to be uploaded
        :param username: NERSC username, needed if different from local username
        :param target_dir: Location where the file should be placed if different from the user default location
        :param session: The requests session to be used or None if the user should be asked to
                         login as part of this call
        :param persist_session: If True then the session will remain open, otherwise we'll log the user
                        out of NERSC at the end
        :param machine: The name of the machine at nersc to be used. Default is 'cori'
        :param register: Boolean indicating whether we should register the file with the openmsi server (default=False)
        :return:
        """
        import requests
        import getpass
        import urllib

        newt_base_url = "https://newt.nersc.gov/newt"
        newt_auth_url = newt_base_url+"/auth"
        newt_logout_url = newt_base_url+"/logout"

        uname = username
        if uname is None:
            print "Please enter the user name"
            if uname == "" or uname is None:
                uname = getpass.getuser()
                tmp = getpass.getpass("Username: (default="+uname+")")
                uname = uname if tmp == "" else tmp
                print tmp, uname
        target = target_dir if target_dir is not None else ("/project/projectdirs/openmsi/omsi_data_private/" + uname)

        s = session
        newt_sessionid = None
        if s is None:  # Login to newt

            s = requests.Session()
            print "Please enter your NERSC Password"
            password = getpass.getpass(prompt="Enter password for user \"" + uname + "\" \n")

            r = s.post(newt_auth_url, {"username": uname, "password": password})
            if r.status_code != 200:
                raise ValueError("Authentication failed.")
            result = r.json()
            if not result['auth']:
                raise ValueError("Authentication failed.")
            newt_sessionid = result['newt_sessionid']

        # Upload the file to NEWT
        upload_url = newt_base_url + "/file/" + machine + target
        fname = os.path.basename(filepath)
        cookies = {"newt_sessionid": newt_sessionid}
        payload = {"name": fname,
                   "file": open(filepath, 'rb'),
                   "submit": "upload"}
        result1 = requests.post(upload_url, files=payload, cookies=cookies)
        if result1.status_code != 200:
            raise ValueError("Upload failed")
        # result = r.json()

        # Set apache ACL at NERSC
        fullname = target + "/" + fname
        acl_url = "https://newt.nersc.gov/newt/command/" + machine
        command = "setfacl -R -m u:48:rwx " + '"' + fullname + '"'
        payload = {"executable": command, "loginenv": "true"}
        result2 = requests.post(acl_url, data=payload, cookies=cookies)

        # Register the file with the openmsi database
        if register:
            add_file_url = os.path.join(WebHelper.default_db_server_url, "openmsi/resources/addfile")
            # addfilepath = fullname.lstrip("/global")
            query_params = {'file': fullname, 'owner': uname}
            add_file_url += "?"
            add_file_url += urllib.urlencode(query_params)
            result_3 = requests.get(url=add_file_url)
        else:
            result_3 = None

        # Logout user if requested
        if not persist_session:
            s = requests.Session()
            r = s.get(newt_logout_url)
            s = None

        return result1, result2, result_3, s

    @staticmethod
    def update_job_status(filepath, db_server, jobid, status='complete'):
        """
        Function used to update the status of the job on the server

        :param filepath: Path of the file to be added to the database (only needed update file permissions)
        :param db_server: The database server url
        :param jobid: The id of the current job.
        :param status: One of 'running', 'complete' or 'error'
        """
        import urllib2
        import urllib

        # If we are at NERSC then set the NERSC Apache permissions
        if 'nersc.gov' in db_server:
            WebHelper.set_apache_acl(filepath)

        # Construct the db add-file url
        update_status_url = os.path.join(db_server, "openmsi/processing/update")
        query_params = {'jobid': jobid, 'status': status}
        update_status_url += "?"
        update_status_url += urllib.urlencode(query_params)

        # Make the url request
        try:
            log_helper.info(__name__, "Updating job status: " + update_status_url)
            url_response = urllib2.urlopen(url=update_status_url)
            if url_response.code == 200:
                return True
        except urllib2.HTTPError as request_error:
            raise ValueError("ERROR: job status could not be updated: \n" +
                             "      Error-code:" + str(request_error.code) + "\n" +
                             "      Error info:" + str(request_error.read()))
        except urllib2.URLError as request_error:
            if sys.version_info >= (2, 7, 9):
                import ssl
                ssl_context = ssl.create_default_context()
                ssl_context.check_hostname = False
                ssl_context.verify_mode = ssl.CERT_NONE
                url_response = urllib2.urlopen(url=update_status_url, context=ssl_context)
                if url_response.code == 200:
                    return True
            else:
                raise ValueError("ERROR: job status could not be updated: \n" +
                                 "      Error-code:" + str(request_error.code) + "\n" +
                                 "      Error info:" + str(request_error.read()))

    @staticmethod
    def register_file_with_db(filepath, db_server, file_user_name, jobid=None, check_add_nersc=True):
        """ Function used to register a given file with the database

            :param filepath: Path of the file to be added to the database
            :param db_server: The database server url
            :param file_user_name: The user to be used, or None if the user should
                                    be determined based on the file URL.
            :param jobid: Optional input parameter defining the jobid to be updated.
                          If the jobid is given then the job will be updated with the
                          database instead of adding the file explicitly. I.e.,
                          instead of register_filer_with_db the update_job_status call
                          is executed.
            :param check_add_nersc: Boolean if set to True performs additional actions to add the
                         file to the pupblic OpenMSI gateway hosted at NERSC.

            :returns: Boolean indicating whether the operation was successful

        """
        import urllib2
        import urllib

        if jobid is not None:
            return WebHelper.update_job_status(filepath=filepath,
                                               db_server=db_server,
                                               jobid=jobid,
                                               status='complete')
        # Check if the
        if db_server == WebHelper.default_db_server_url and check_add_nersc:
            is_allowed_path = False
            for allowed_nersc_location in WebHelper.allowed_nersc_locations:
                if filepath.startswith(allowed_nersc_location):
                    is_allowed_path = True
                    break
            if not is_allowed_path and file_user_name in WebHelper.super_users:
                print "WARNING: Attempt to add a file to openmsi.nersc.gov that is not in a default location."
                print "Do you want to add the file? (Y/N):"
                num_trys = 3
                timeout = 5*60  # Timeout after 5 minutes
                for i in range(num_trys):
                    # user_input = raw_input()
                    user_input = UserInput.userinput_with_timeout(timeout=timeout, default=None)
                    if user_input is None:
                        warnings.warn("WARNING: Attempt to add a file to openmsi.nersc.gov that," +
                                      " is not in a default location. Timeout occurred before" +
                                      " user confirmed. Aborted adding the file to the DB.")
                        return False
                    if user_input == "Y" or user_input == "y" or user_input == "Yes" or \
                            user_input == "yes" or user_input == "YES":
                        break
                    elif user_input == "N" or user_input == "n" or user_input == "No" or \
                            user_input == "no" or user_input == "NO":
                        return False
                    else:
                        if i == (num_trys - 1):
                            warnings.warn("WARNING: Attempt to add a file to openmsi.nersc.gov that," +
                                          " is not in a default location. User input unrecognized." +
                                          " Aborted adding the file to the DB.")
                            return False
                        print "Unrecognized response. Do you want to add the file? (Y/N): "
            elif not is_allowed_path:
                warnings.warn("Adding file to the OpenMSI database in unconventional location not permitted for user.")
                return False
            else:
                pass  # Adding the file to the db is allowed

        # If we are at NERSC then set the NERSC Apache permissions
        if 'nersc.gov' in db_server and check_add_nersc:
            WebHelper.set_apache_acl(filepath)

        # Determine the user
        curr_user = file_user_name
        if not curr_user:
            curr_user = os.path.dirname(filepath).split("/")[-1]
        if not curr_user:
            raise ValueError("ERROR: File could not be added to DB. Owner could not be determined.")

        # Construct the db add-file url
        add_file_url = os.path.join(db_server, "openmsi/resources/addfile")
        addfilepath = filepath
        # Correct the filepath if we are on openmsi.nersc.gov, as /global is not mounted but only /project.
        if db_server == WebHelper.default_db_server_url and addfilepath.startswith("/global/project/projectdirs"):
            addfilepath = filepath.lstrip("/global")
        query_params = {'file': os.path.abspath(addfilepath), 'owner': curr_user}
        add_file_url += "?"
        add_file_url += urllib.urlencode(query_params)
        # add_file_url = add_file_url + "?file=" + \
        #    os.path.abspath(filepath) + "&user=" + curr_user

        # Make the url request
        try:
            log_helper.info(__name__, "Registering file with DB: " + add_file_url)
            url_response = urllib2.urlopen(url=add_file_url)
            if url_response.code == 200:
                return True
        except urllib2.HTTPError as request_error:
            raise ValueError("ERROR: File could not be added to DB: \n" +
                             "      Error-code:" + str(request_error.code) + "\n" +
                             "      Error info:" + str(request_error.read()))

        return False

    @staticmethod
    def set_apache_acl(filepath):
        """
        Helper function used to set acl permissions for apache to make the given file accesible
        to Apache at NERSC. This necessary to make the file readable for adding it to the
        database.

        :param filepath: String with the path to the file for which ACL permission should be set

        """
        log_helper.info(__name__, "Setting NERSC ACL permissions for Apache")
        # Note u:48 is a replacement for u:apache to ensure that
        # that the command works properly on edison.nersc.gov which
        # does not have the apache user. However u:48 is equivalent.
        command = "setfacl -R -m u:48:rwx " + '"' + filepath + '"'
        os.system(command)

    @staticmethod
    def send_email(subject,
                   body,
                   sender='convert@openmsi.nersc.gov',
                   email_type='success',
                   email_success_recipients=None,
                   email_error_recipients=None):
        """
        Send email notification to users.

        :param subject: Subject line of the email
        :param body: Body text of the email.
        :param sender: The originating email address
        :param email_type: One of 'success, 'error', 'warning'. Error messages are sent
                     to ConvertSettings.email_error_recipients, success messages to
                     ConvertSettings.email_success_recipients and warning messages are sent to both lists.
        :param email_success_recipients: List of user that should receive an email if the status is success
                     or warning.
        :param email_error_recipients: List of users that should receive an email if the status is error
                     or warning.

        """
        if email_error_recipients is None:
            email_error_recipients = []
        if email_success_recipients is None:
            email_success_recipients = []
        # Define the list of recipients
        if email_type == 'success':
            recipients = email_success_recipients
        elif email_type == 'error':
            recipients = email_error_recipients
        else:
            recipients = email_error_recipients + email_success_recipients
        # Remove duplicates from the list of recipients
        recipients = list(set(recipients))
        # Check if we have any recipients
        if len(recipients) == 0:
            return

        from smtplib import SMTP
        from email.MIMEText import MIMEText
        from email.Header import Header
        from email.Utils import parseaddr, formataddr

        header_charset = 'ISO-8859-1'
        body_charset = 'US-ASCII'
        for current_char_set in 'US-ASCII', 'ISO-8859-1', 'UTF-8':
            try:
                body_charset = current_char_set
                body.encode(body_charset)
            except UnicodeError:
                pass
            else:
                break

        # Define the sender and recipients
        sender_name, sender_addr = parseaddr(sender)
        sender_name = str(Header(unicode(sender_name), header_charset))
        sender_addr = sender_addr.encode('ascii')

        tostr = ""
        for ri in range(len(recipients)):
            rec = recipients[ri]
            recname, recaddr = parseaddr(rec)
            recname = str(Header(unicode(recname), header_charset))
            recaddr = recaddr.encode('ascii')
            tostr += formataddr((recname, recaddr))
            if ri < (len(recipients)-1):
                tostr += ", "

        # Construct the message
        msg = MIMEText(body.encode(body_charset), 'plain', body_charset)
        msg['From'] = formataddr((sender_name, sender_addr))
        msg['To'] = tostr
        msg['Subject'] = Header(unicode(subject), header_charset)

        # Send the message using sendmail
        try:
            smtp = SMTP("localhost")
            smtp.sendmail(sender, recipients, msg.as_string())
            smtp.quit()
        except:
            warnings.warn('Email could not be sent' + str(sys.exc_info()))

    """
    def loginUser(requestPassword=False):

        import getpass
        import urllib2
        import urllib
        #Setup the user login mechanism in urllib2
        if ConvertSettings.db_server_url.startswith("http:"):
            ConvertSettings.db_server_url = "https:" + ConvertSettings.db_server_url.lstrip("http:")
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        username = raw_input("Login: >> ") #Enter username
        password = getpass.getpass()       #Enter password
        password_mgr.add_password(None, ConvertSettings.db_server_url, username, password)
        handler = urllib2.HTTPBasicAuthHandler(password_mgr)
        opener = urllib2.build_opener(handler)
        urllib2.install_opener(opener)

        #Login the user to the website
        login_data = urllib.urlencode({
                'username' : username,
                'password' : password,
            })

        response = urllib2.urlopen(ConvertSettings.db_server_url, login_data)
        print response
    """


####################################################################
####################################################################
# User input helper functions                                      #
####################################################################
####################################################################
class UserInput(object):
    """Collection of helper functions used to collect user input"""
    def __init__(self):
        pass

    @staticmethod
    def userinput_with_timeout_default(timeout, default=''):
        """Read user input. Return default value given after timeout.

          :param timeout: Number of seconds till timeout
          :param default: Default string to be returned after timeout
          :type default: String

          :returns: String

        """
        sys.stdout.flush()
        rlist, _, _ = select([sys.stdin], [], [], timeout)
        if rlist:
            userinput = sys.stdin.readline().replace('\n', '')
        else:
            userinput = default
        return userinput

    @staticmethod
    def userinput_with_timeout_windows(timeout, default=''):
        """Read user input. Return default value given after timeout.
           This function is used when running on windows-based systems.

          :param timeout: Number of seconds till timeout
          :param default: Default string to be returned after timeout
          :type default: String

          :returns: String

        """
        start_time = time.time()
        sys.stdout.flush()
        userinput = ''
        while True:
            if msvcrt.kbhit():
                readchar = msvcrt.getche()
                if ord(readchar) == 13:  # enter_key
                    break
                elif ord(readchar) >= 32:  # space_char
                    userinput += readchar
            if len(userinput) == 0 and (time.time() - start_time) > timeout:
                break
        if len(userinput) > 0:
            return userinput
        else:
            return default

    @staticmethod
    def userinput_with_timeout(timeout, default=''):
        """
        Read user input. Return default value given after timeout.
        This function decides which platform-dependent version should
        be used to retrieve the user input.

        :param timeout: Number of seconds till timeout
        :param default: Default string to be returned after timeout
        :type default: String

        :returns: String
        """
        if platform.system() == "Windows":
            return UserInput.userinput_with_timeout_windows(timeout, default)
        else:
            return UserInput.userinput_with_timeout_default(timeout, default)
