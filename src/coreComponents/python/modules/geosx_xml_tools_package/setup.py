from distutils.core import setup

setup(
    name="geosx_xml_tools",
    version="0.5.0",
    description="Tools for enabling advanced xml features in GEOSX",
    author="Chris Sherman",
    author_email="sherman27@llnl.gov",
    packages=["geosx_xml_tools", "geosx_xml_tools.tests"],
    entry_points={
        "console_scripts": [
            "preprocess_xml = geosx_xml_tools.main:preprocess_serial",
            "format_xml = geosx_xml_tools.xml_formatter:main",
            "test_geosx_xml_tools = geosx_xml_tools.tests.test_manager:main",
            "test_attribute_coverage = geosx_xml_tools.attribute_coverage:main",
            "test_xml_redundancy = geosx_xml_tools.xml_redundancy_check:main",
        ]
    },
    install_requires=["lxml>=4.5.0", "parameterized"],
)
