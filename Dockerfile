FROM python:3.11-ubi8

# 1. Start as root to perform system-level changes
USER root

# 2. Install the spatial library (Confirmed working name)
RUN yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    yum install -y spatialindex spatialindex-devel && \
    yum clean all

# 3. Copy your code into the temporary staging area while still root
COPY . /tmp/src

# 4. NOW that the files exist, give ownership to the OpenShift user (1001)
# We also give ownership of the destination folder (/opt/app-root/src)
RUN chown -R 1001:0 /tmp/src /opt/app-root/src

# 5. Switch to the unprivileged user for safety and S2I compatibility
USER 1001

# 6. Run the assembly script. Since 1001 owns /tmp/src, it can now move files.
RUN /usr/libexec/s2i/assemble

# 7. Set the startup command
CMD /usr/libexec/s2i/run