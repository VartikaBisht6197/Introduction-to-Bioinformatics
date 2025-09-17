FROM mambaorg/micromamba:1.5.1

# Copy environment file
COPY Intro_to_Bioinformatics.yml /tmp/Intro_to_Bioinformatics.yml

# Install everything into base env
RUN micromamba install -y -f /tmp/Intro_to_Bioinformatics.yml -n base && \
    micromamba clean --all -f -y

# Initialize micromamba for all future shells
SHELL ["bash", "-c"]
RUN micromamba shell init -s bash -p /opt/conda

