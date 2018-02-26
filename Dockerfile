FROM amirkdv/science-py-env

RUN apt-get update -qq

# Dependencies for python-cffi-backend
RUN apt-get install -qqy libffi6 libffi-dev

# Dependencies for sqlite3 and the CLI for debugging
RUN apt-get install -qqy libsqlite3-dev sqlite3

# install doxygen and graphviz for docs
RUN apt-get install -qqy doxygen graphviz

# install HTS for pysam
RUN apt-get install -qqy libhts1
