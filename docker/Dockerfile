FROM ubuntu:20.04

RUN apt-get update -y && \
  apt-get install -y gcc libblas-dev liblapacke-dev make vim

WORKDIR /app
COPY ../src /app/src
COPY ../include /app/include
COPY ../Makefile /app/
COPY ../docker.mk /app/

ENV HOSTNAME=docker
RUN mkdir -p bin && \
    ls -la && \
    make
