import sim_launcher
from os import environ

# cf. https://stackoverflow.com/questions/48066106/converting-vector-looking-string-into-vector-values-python
lmin=[int(i.strip()) for i in environ.get('ADAPTATION_MINIMUM_LEVEL')[1:-1].split(",")]

yIsX = environ.get('ADAPTATION_PARAM_Y_EQUALS_X', 'False').lower() in ['true', '1']
if yIsX:
    assert(lmin[0]==lmin[1])
print("Starting simulation at lowest resolution " + str(lmin))
sim_launcher.dispatch_to_run(level_vector=lmin)

# try starting a whole bunch more, to save on overall waiting time
# at least start the direct neighbors
for d in range(len(lmin)):
    neighbor = lmin
    neighbor[d] += 1
    print("Starting simulation at resolution " + str(neighbor))
    sim_launcher.dispatch_to_run(level_vector=neighbor)

# then assume that x is refined at least twice
for i in [1,2]:
    x_high = lmin
    x_high[0] += i
    for d in range(len(x_high)):
        neighbor = x_high
        neighbor[d] += 1
        try:
            sim_launcher.dispatch_to_run(level_vector=neighbor)
            print("Starting simulation at resolution " + str(neighbor))
        except sim_launcher.ProbFolderExistsError as e:
            pass
# if y is not x, assume it is refined at least twice, too
if not yIsX:
    for i in [1,2]:
        y_high = lmin
        y_high[1] += i
        for d in range(len(y_high)):
            neighbor = y_high
            neighbor[d] += 1
            try:
                sim_launcher.dispatch_to_run(level_vector=neighbor)
                print("Starting simulation at resolution " + str(neighbor))
            except sim_launcher.ProbFolderExistsError as e:
                pass

# and z is refined at least once
z_high = lmin
z_high[2] += 1
for d in range(len(z_high)):
    neighbor = z_high
    neighbor[d] += 1
    try:
        sim_launcher.dispatch_to_run(level_vector=neighbor)
        print("Starting simulation at resolution " + str(neighbor))
    except sim_launcher.ProbFolderExistsError as e:
        pass