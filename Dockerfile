# Use an official Ubuntu as a parent image
FROM ubuntu:20.04

# Set environment variables to avoid user interaction during installation
ENV DEBIAN_FRONTEND=noninteractive

# Update and install dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \
    libgsl0-dev \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

# Clone the npstat repository
RUN git clone https://github.com/lucaferretti/npstat.git . 

# Run the autoreconf, configure, make, and make install commands
RUN make 
WORKDIR /data
# Set the entrypoint to run the npstat binary
ENTRYPOINT ["/app/npstat"]

# By default, run npstat with no arguments
CMD ["--help"]
