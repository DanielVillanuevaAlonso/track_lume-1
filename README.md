# track_lume-1
This web app tracks and monitors the LUME-1 satellite, automating the calculation of future satellite passes and connecting to the antenna rotor controller in a ground station.

# Deployment
First, build the Docker image:

	docker build -t <tag_name> .

Then, deploy the Docker container:

	docker run -p 5000:5000 <tag_name>
