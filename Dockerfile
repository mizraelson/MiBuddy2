FROM ubuntu

ARG mixcrVersion=3.0.12
ARG migecVersion=1.2.9
ARG vdjtoolsVersion=1.2.1

RUN apt-get update
RUN apt-get install -y wget python3-pip unzip default-jdk
RUN rm -r -f /var/lib/apt/lists/*
RUN mkdir /data
RUN mkdir /MiBuddy2

RUN cd / \
    && wget https://github.com/milaboratory/mixcr/releases/download/v${mixcrVersion}/mixcr-${mixcrVersion}.zip \
    && unzip mixcr-${mixcrVersion}.zip \
    && mv mixcr-${mixcrVersion} mixcr \
    && rm mixcr-${mixcrVersion}.zip

RUN mkdir /migec \
    && cd migec/ \
    && wget https://github.com/mikessh/migec/releases/download/${migecVersion}/migec-${migecVersion}.zip \
    && unzip migec-${migecVersion}.zip \
    && rm migec-${migecVersion}.zip \
    && rm migec


RUN cd / \
    && wget https://github.com/mikessh/vdjtools/releases/download/${vdjtoolsVersion}/vdjtools-${vdjtoolsVersion}.zip \
    && unzip vdjtools-${vdjtoolsVersion}.zip \
    && mv vdjtools-${vdjtoolsVersion} vdjtools \
    && rm vdjtools-${vdjtoolsVersion}.zip \
    && rm /vdjtools/vdjtools

COPY MiBuddy2.py requirements.txt /MiBuddy2/
COPY migec /migec/
COPY vdjtools /vdjtools/

RUN pip3 install -r /MiBuddy2/requirements.txt \
    && chmod +x /MiBuddy2/MiBuddy2.py \
    && chmod +x /vdjtools/vdjtools \
    && chmod +x /migec/migec

ENV PATH="/mixcr:/migec:/vdjtools:/MiBuddy2:${PATH}"
ENV JAVA_XMX="7g"

WORKDIR /data
ENTRYPOINT ["MiBuddy2.py"]


