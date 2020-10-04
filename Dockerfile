FROM python:2.7-slim

# for working inside container
RUN apt update -y && apt install vim -y

WORKDIR /app

RUN pip install scipy numpy

COPY . /app

RUN pip install -e .

CMD ["python", "examples/psychotria.lagrange.py"]