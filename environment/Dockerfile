# hash:sha256:9cc8eee93f808b469c58026d07522b116d361fcad44a95b19d96f018e248fcc2
FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && add-apt-repository -y 'deb http://deb.codeocean.com/rstudio-server-bionic/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common
