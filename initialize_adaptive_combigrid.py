import sim_launcher
from os import environ

# cf. https://stackoverflow.com/questions/48066106/converting-vector-looking-string-into-vector-values-python
lmin=[int(i.strip()) for i in environ.get('ADAPTATION_MINIMUM_LEVEL')[1:-1].split(",")]

yIsX = os.environ.get('ADAPTATION_PARAM_Y_EQUALS_X', 'False').lower() in ['true', '1']
if yIsX:
    assert(lmin[0]==lmin[1])
print("Starting simulation at lowest resolution " + str(lmin))
sim_launcher.dispatch_to_run(level_vector=lmin)