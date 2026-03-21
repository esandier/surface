FROM python:3.11-ubi8

USER root
# Install the library
RUN yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    yum install -y spatialindex spatialindex-devel && \
    yum clean all

# Ditch the S2I 'assemble' script—it's too picky about permissions.
# We will do a standard Python install ourselves.
WORKDIR /opt/app-root/src
COPY . .
RUN chown -R 1001:0 /opt/app-root/src
USER 1001
RUN pip install --no-cache-dir -r requirements.txt

CMD python manage.py runserver 0.0.0.0:8080