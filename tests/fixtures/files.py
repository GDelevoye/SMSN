import pytest
import smsn
import os


@pytest.fixture
def get_resource():
    """Input = A relative path from the resource directory
    Output = An absolute path to the resource in question

    """
    resources_dir = smsn.pipeline._getAbsPath("/resources/")
    return resources_dir