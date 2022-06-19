FROM ubuntu:20.04

MAINTAINER Daniel Villanueva Alonso

EXPOSE 5000
RUN apt update -y \ 
    && apt upgrade -y \
	&& apt install -y python3-pip python3-dev

COPY requirements.txt .
COPY  AppWebTrackLume-1/ AppWebTrackLume-1/
RUN pip3 install -r requirements.txt
WORKDIR "AppWebTrackLume-1"
CMD python3 main.py runserver 0.0.0.0:5000