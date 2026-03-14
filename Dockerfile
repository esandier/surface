FROM python:3.11-ubi8

USER root
# 1. Install the system library
RUN yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    yum install -y spatialindex spatialindex-devel && \
    yum clean all

# 2. Fix permissions so the S2I user (1001) can move files
RUN chown -R 1001:0 /opt/app-root/src /tmp/src

# 3. Switch back to the standard unprivileged user
USER 1001

# 4. Standard S2I assembly
COPY . /tmp/src
RUN /usr/libexec/s2i/assemble

# 5. Set the default command
CMD /usr/libexec/s2i/run
