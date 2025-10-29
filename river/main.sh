cd analysis

# Download the database files
# wget https://zenodo.org/records/10394426/files/CircNetVis_resources.tar.gz?download=1 -O CircNetVis_resources.tar.gz
# tar -xvzf CircNetVis_resources.tar.gz -C CircNetVis --strip-components 1
if [ -d "CircNetVis_mm9/Tools" ]; then
    echo "CircNetVis_mm9 resources already exist."
else
    wget https://zenodo.org/records/10394426/files/CircNetVis_mm9_resources.tar.gz?download=1 -O CircNetVis_mm9_resources.tar.gz
    tar -xvzf CircNetVis_mm9_resources.tar.gz -C CircNetVis_mm9 --strip-components 1
fi
# Run the Shiny app
pixi run Rscript -e "shiny::runApp('CircNetVis_mm9')"