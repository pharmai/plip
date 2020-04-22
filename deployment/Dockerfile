FROM ubuntu:18.04 AS builder

ARG OPENBABEL_TAG="openbabel-3-0-0"
ARG RDKIT_TAG="Release_2019_09_3"
ARG MMTF_TAG="v1.0.0"
ARG PYMOL_TAG="v2.3.0"

LABEL maintainer="PharmAI GmbH <contact@pharm.ai>" \
        org.label-schema.name="PLIP: The Protein-Ligand Interaction Profiler" \
        org.label-schema.description="https://www.doi.org/10.1093/nar/gkv315"

RUN apt-get update && apt-get install -y \
    cmake \
    git \
    g++ \
    imagemagick \
    libboost-all-dev \
    libeigen3-dev \
    libfreetype6-dev \
    libghc-zlib-dev \
    libglew-dev \
    libglm-dev \
    libmsgpack-dev \
    libnetcdf-dev \
    libxml2-dev \
    libpng-dev \
    libpython3-all-dev \
    python3 \
    python3-numpy \
    python3-pyqt5.qtopengl \
    swig

# build OpenBabel
WORKDIR /src
RUN git clone -b $OPENBABEL_TAG \
    https://github.com/openbabel/openbabel.git
WORKDIR /src/openbabel/build
RUN cmake .. \
    -DPYTHON_EXECUTABLE=/usr/bin/python3.6 \
    -DPYTHON_BINDINGS=ON \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DRUN_SWIG=ON
RUN make -j$(nproc --all) install

# build RDkit
WORKDIR /src
RUN git clone -b $RDKIT_TAG \
    https://github.com/rdkit/rdkit.git
WORKDIR /src/rdkit/build
RUN cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DRDK_BUILD_INCHI_SUPPORT=ON \
    -DPYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.6m.so \
    -DPYTHON_INCLUDE_DIR=/usr/include/python3.6 \
    -DPYTHON_EXECUTABLE=/usr/bin/python3.6
RUN make -j$(nproc --all) install
ENV RDBASE /src/rdkit/
ENV LD_LIBRARY_PATH /usr/local/lib/:$RDBASE/lib
ENV PYTHONPATH $PYTHONPATH:$RDBASE

# build mmtf-cpp
WORKDIR /src
RUN git clone -b $MMTF_TAG \
    https://github.com/rcsb/mmtf-cpp.git
WORKDIR /src/mmtf-cpp/build
RUN cmake ..
RUN make -j$(nproc --all) install

# build PyMOL
WORKDIR /src
RUN git clone -b $PYMOL_TAG \
    https://github.com/schrodinger/pymol-open-source
WORKDIR /src/pymol-open-source
RUN python3 setup.py build install \
    --prefix=/usr/local \
    --jobs=$(nproc --all)