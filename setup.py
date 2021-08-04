import setuptools

try:
    from pyfiglet import Figlet
    f = Figlet(font='speed')
    long_message = f.renderText('TauRunner')
except:
    long_message = 'TauRunner'

version = "0.0.4"

setuptools.setup(
    name="taurunner", 
    version=version,
    author="Safa, I. et al.",
    author_email="isafa@wisc.edu",
    description="Code for propagating ultra-high-energy neutrinos",
    long_description=long_message,
    #long_description_content_type="text/markdown",
    url="https://github.com/icecube/TauRunner",
    packages=setuptools.find_packages(),
    install_requires=['numpy>=1.16.6',
                      'scipy>=1.2.3',
                      'proposal==6.1.6'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.1',
    package_data={"taurunner.resources.solar_models"         : ["*.txt"],
                  "taurunner.resources.cross_section_tables" : ["*"],
                  "taurunner.resources.secondaries_splines"  : ["*"],
                  "taurunner.resources.proposal_tables"      : ["tables.txt"],}    
)
