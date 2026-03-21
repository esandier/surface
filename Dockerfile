FROM image-registry.openshift-image-registry.svc:5000/openshift/python:3.11

USER root
# Installation des dépendances système nécessaires pour Rtree
RUN yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    yum install -y libspatialindex libspatialindex-devel && \
    yum clean all

WORKDIR /opt/app-root/src
COPY . .
RUN chown -R 1001:0 /opt/app-root/src

USER 1001
RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "manage.py", "runserver", "0.0.0.0:8080"]