import glob, os

listOfFiles = []

for (dirpath, dirnames, filenames) in os.walk('./resultado/'):
    listOfFiles += [os.path.join(dirpath, file) for file in filenames]

for elem in listOfFiles:
    if elem[-4:] == 'gpkg':
        os.rename(elem, elem.replace(' - ', '_').replace(' ', '_'))
        # print(elem.replace(' - ', '_').replace(' ', '_'))