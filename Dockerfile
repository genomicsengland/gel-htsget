FROM python:2.7
RUN apt-get update && apt-get install -y \
    build-essential \
    ntp \
	libpq-dev \
	python-dev \
	git \
	supervisor \
	nginx

ENV PYTHONUNBUFFERED 1
RUN mkdir -p /app/htsget_server
WORKDIR /app
ADD ./htsget_server /app/htsget_server/
ADD ./setup.py ./VERSION /app/
RUN pip install . -i https://pypi.gel.zone/genomics/dev

COPY ./config/supervisord.conf /etc/supervisor/supervisord.conf
COPY ./config/opencga.json /etc/htsget_server/opencga.json
COPY ./config/vcf_types.tsv /etc/vcf_types.tsv

CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/supervisord.conf", "-n"]