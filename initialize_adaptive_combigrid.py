import sim_launcher
from os import environ

# cf. https://stackoverflow.com/questions/48066106/converting-vector-looking-string-into-vector-values-python
lmin=[int(i.strip()) for i in environ.get('ADAPTATION_MINIMUM_LEVEL')[1:-1].split(",")]

print("Starting simulation at lowest resolution " + str(lmin))
sim_launcher.dispatch_to_run(level_vector=lmin)