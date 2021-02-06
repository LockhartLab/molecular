
with open('requirements.txt', 'r') as buffer:
    data = buffer.readlines()

with open('pyproject.toml', 'w') as buffer:
    buffer.write('[build-system]\n')
    buffer.write('# note: automatically generated\n')
    buffer.write('requires = [\n')
    for line in data:
        buffer.write('    \'line\'\n')
    buffer.write(']')