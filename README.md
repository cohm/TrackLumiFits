# TrackLumiFits

Boilerplate for toy simulation studies testing how well the luminosity measurement using track counting can be improved by fitting a function describing the number of tracks per bunch crossing.

# How to use this
This repository contains a starting point for how to use the `ROOT` library from `python` in a jupyter notebook. There are two main options for doing this:
* install jupyter notebook directly on your own computer
* use it from a docker container

For the latter option, install docker and put this in your terminal:

```
docker run -v /Users/cohm/Desktop/tmp/TrackLumiFits/:/work -p 3000:8080 cohm/pyroot-notebook:latest
```

This will start the docker image `cohm/pyroot-notebook` which contains a jupyter notebook server, and make the path after the `-v` switch (which could e.g. be the path to where you checked out this repository) accessible from inside the notebooks. When starting this image you should get output like this in your terminal:

```
OhmBookPro13:~/Desktop/tmp > docker run -v /Users/cohm/Desktop/tmp/TrackLumiFits/:/work -p 3000:8080 cohm/pyroot-notebook:latest
(...)
[I 21:07:22.848 NotebookApp] Writing notebook server cookie secret to /root/.local/share/jupyter/runtime/notebook_cookie_secret
[I 21:07:23.049 NotebookApp] Serving notebooks from local directory: /work
[I 21:07:23.050 NotebookApp] Jupyter Notebook 6.1.4 is running at:
[I 21:07:23.051 NotebookApp] http://f4a437845500:8080/?token=2b9ef203e4805fce2baf10e4022a01388336c0c90ddec870
[I 21:07:23.051 NotebookApp]  or http://127.0.0.1:8080/?token=2b9ef203e4805fce2baf10e4022a01388336c0c90ddec870
[I 21:07:23.051 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[W 21:07:23.055 NotebookApp] No web browser found: could not locate runnable browser.
[C 21:07:23.055 NotebookApp] 
    
    To access the notebook, open this file in a browser:
        file:///root/.local/share/jupyter/runtime/nbserver-1-open.html
    Or copy and paste one of these URLs:
        http://f4a437845500:8080/?token=2b9ef203e4805fce2baf10e4022a01388336c0c90ddec870
     or http://127.0.0.1:8080/?token=2b9ef203e4805fce2baf10e4022a01388336c0c90ddec870

```

If you then open the link `http://127.0.0.1:3000/?token=2b9ef203e4805fce2baf10e4022a01388336c0c90ddec870` in a browser on your laptop, you can start coding. 

(NB! The port number 8080 in the terminal print-out is changed to 3000 here)
