# Start from the same Python image you are currently using
FROM registry.access.redhat.com/ubi8/python-311

USER root

# 1. Install the EPEL repository (where libspatialindex lives)
RUN yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm && \
    # Disable all repos, then only enable the ones we need to save RAM
    yum install -y --disablerepo="*" --enablerepo="ubi-8-appstream-rpms" --enablerepo="ubi-8-baseos-rpms" --enablerepo="epel" libspatialindex && \
    yum clean all

# 2. Switch back to the standard OpenShift user
USER 1001

# 3. Trigger the standard S2I assemble script to install your requirements.txt
COPY . /tmp/src
RUN /usr/libexec/s2i/assemble

# 4. Set the default startup command
CMD /usr/libexec/s2i/run