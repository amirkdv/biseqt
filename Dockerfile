FROM amirkdv/science-py-env

RUN apt-get update -qq

# Dependencies for python-cffi-backend
RUN apt-get install -qqy libc6 libffi6 libffi-dev

# Dependencies for python-igraph
RUN apt-get install -qqy libc6 libigraph0v5 libigraph0-dev

# libcairo and python cairo bindings are needed for igraph plotting
# NOTE pycairo is no good, instead use cairocffi.
RUN apt-get install -qqy libcairo2

# Clear the font cache
RUN fc-cache -vf

# install doxygen and pandoc for docs
RUN apt-get install -qqy doxygen pandoc
