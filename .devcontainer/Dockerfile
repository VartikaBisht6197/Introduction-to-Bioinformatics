FROM mambaorg/micromamba:1.5.1

# Copy environment file
COPY Intro_to_Bioinformatics.yml /tmp/Intro_to_Bioinformatics.yml

# Install everything into base env
RUN micromamba create -y -n main -f /tmp/Intro_to_Bioinformatics.yml

# Initialize micromamba for all future shells
SHELL ["bash", "-c"]
RUN micromamba shell init -s bash -p /opt/conda

