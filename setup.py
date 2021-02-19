import setuptools

try:
    from pyfiglet import Figlet
    f = Figlet(font='speed')
    long_message = f.renderText('TauRunner')
except:
    long_message = 'TauRunner'

version = "0.0.1"

setuptools.setup(
    name="taurunner", # Replace with your own username
    version=version,
    author="Safa, I.",
    author_email="isafa@wisc.edu",
    description="Code for propagating ultra-high-energy neutrinos",
    long_description=long_message,
    #long_description_content_type="text/markdown",
    url="https://github.com/icecube/TauRunner",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.1',
)