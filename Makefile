build:
	docker build -t lagrange .

run:
	docker run -it lagrange bash

test:
	docker run -it lagrange