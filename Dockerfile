# Start from the same Python image you are currently using
FROM python:3.11-ubi8
USER root

# 1. Install the EPEL repository (where libspatialindex lives)
USER root
RUN yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    yum install -y 'dnf-command(config-manager)' && \
    yum config-manager --set-enabled epel && \
    yum install -y libspatialindex-devel && \
    yum clean all

# 2. Switch back to the standard OpenShift user
USER 1001

# 3. Trigger the standard S2I assemble script to install your requirements.txt
COPY . /tmp/src
RUN /usr/libexec/s2i/assemble

# 4. Set the default startup command
CMD /usr/libexec/s2i/run