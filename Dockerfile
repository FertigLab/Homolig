#add --platform=linux/amd64 to docker build command to build on a Mac M1
FROM python:3.11

COPY . /homolig
WORKDIR /homolig

RUN apt-get update && apt-get install build-essential patch -y
RUN pip install --upgrade pip 

RUN python ./setup.py install
