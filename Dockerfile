FROM ubuntu:latest

RUN mkdir /code
WORKDIR /code

env DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y build-essential python \
    python-dev python-pip libpq-dev \
    libsasl2-dev libldap2-dev libssl-dev git \
    libz-dev libffi-dev
ENV PYTHONUNBUFFERED 1
ADD . /code

RUN pip install --upgrade pip==9.0.3 && pip install -r /code/requirements.txt --index-url=https://pypi.gel.zone/genomics/dev

CMD [ "python", "/code/htsget_server/htsget_app.py" ]

EXPOSE 8888