FROM mambaorg/micromamba:1.5.1

# Copy environment file
COPY Intro_to_Bioinformatics.yml /tmp/Intro_to_Bioinformatics.yml

# Install everything into base env
RUN micromamba install -y -n main -f /tmp/Intro_to_Bioinformatics.yml

# Make R kernel visible to Jupyter
RUN R -e "IRkernel::installspec(user = FALSE)"
