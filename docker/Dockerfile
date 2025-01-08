FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Paris

RUN apt-get update -y && \
    apt-get install -y gcc libblas-dev liblapacke-dev liblapack-dev make vim gnuplot tzdata && \
    ln -fs /usr/share/zoneinfo/Europe/Paris /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata

WORKDIR /app
COPY . .

ENV HOSTNAME=docker
