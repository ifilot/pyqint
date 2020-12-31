FROM quay.io/pypa/manylinux2010_x86_64
RUN yum -y install eigen3-devel boost-devel
RUN yum -y install mlocate
RUN yum -y install gcc
RUN updatedb
RUN /opt/python/cp35-cp35m/bin/python -m pip install numpy tqdm cython
RUN /opt/python/cp36-cp36m/bin/python -m pip install numpy tqdm cython
RUN /opt/python/cp37-cp37m/bin/python -m pip install numpy tqdm cython
RUN /opt/python/cp38-cp38/bin/python -m pip install numpy tqdm cython
RUN /opt/python/cp39-cp39/bin/python -m pip install numpy tqdm cython
