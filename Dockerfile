ARG DEALII_IMAGE_VERSION="master"

FROM dealii/dealii:${DEALII_IMAGE_VERSION}-noble AS builder

USER root

WORKDIR /opt/adaflo

COPY . .

RUN rm -rf build && \
    cmake \
      -S . \
      -B build \
      -D CMAKE_CXX_STANDARD=20 \
      -D BUILD_SHARED_LIBS=ON \
      -D CMAKE_BUILD_TYPE=Release && \
    cmake --build build -j$(nproc) --target adaflo

ENV ADAFLO_LIB=/opt/adaflo/build/
ENV ADAFLO_INCLUDE=/opt/adaflo/include/

WORKDIR /workspace
