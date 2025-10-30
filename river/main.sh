cd analysis

# Download the database files
if [ -d "CircNetVis/Tools" ]; then
    echo "Database already exists."
else
    wget https://zenodo.org/records/10394426/files/CircNetVis_database.tar.gz?download=1 -O CircNetVis_database.tar.gz
    tar -xvzf CircNetVis_database.tar.gz -C CircNetVis --strip-components 1
fi
# Run the Shiny app
pixi run Rscript -e "shiny::runApp('CircNetVis_mm9')"