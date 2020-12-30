FROM quay.io/pypa/manylinux2014_x86_64
RUN yum -y install eigen3-devel boost-devel
RUN yum -y install mlocate
RUN yum -y install gcc
RUN updatedb
