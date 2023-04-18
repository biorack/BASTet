FROM mambaorg/micromamba:1.4.2
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

WORKDIR /home/mambauser
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN git clone https://github.com/biorack/BASTet.git
ENV PYTHONPATH=/home/mambauser/BASTet

CMD ["/opt/conda/bin/python"]
