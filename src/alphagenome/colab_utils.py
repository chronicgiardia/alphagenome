# Copyright 2025 Google LLC.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Utility functions for Google Colab."""

import os
import pandas as pd

# Check if we're running in Google Colab
try:
  from google.colab import data_table, files, userdata as _userdata
  IN_COLAB = True
  # Enable dataframe formatter if in Colab
  data_table.enable_dataframe_formatter()
  data_table.DataTable.max_rows = 100_000
except ImportError:
  IN_COLAB = False
  # Set pandas options for better local display
  pd.set_option('display.max_rows', 100_000)
  pd.set_option('display.max_columns', None)
  pd.set_option('display.width', None)
  pd.set_option('display.max_colwidth', None)
  
  # Create mock objects for Colab-specific functionality
  class MockDataTable:
    max_rows = 100_000
    
    @staticmethod
    def enable_dataframe_formatter():
      pass
  
  class MockFiles:
    @staticmethod
    def download(filename):
      print(f"File '{filename}' would be downloaded in Google Colab. ")
      print(f"In local environment, you can find the file at: {filename}")
    
    @staticmethod
    def upload():
      print("File upload is only available in Google Colab.")
      print("In local environment, please specify file paths directly.")
      return {}
  
  data_table = MockDataTable()
  files = MockFiles()


def get_api_key(secret: str = 'ALPHA_GENOME_API_KEY'):
  """Returns API key from environment variable or Colab secrets.

  Tries to retrieve the API key from the environment first. If not found,
  attempts to retrieve it from Colab secrets (if running in Colab).

  Args:
    secret: The name of the environment variable or Colab secret key to
      retrieve.

  Raises:
    ValueError: If the API key cannot be found in the environment or Colab
      secrets.
  """

  if api_key := os.environ.get(secret):
    return api_key

  if IN_COLAB:
    try:
      api_key = _userdata.get(secret)
      return api_key
    except (
        _userdata.NotebookAccessError,
        _userdata.SecretNotFoundError,
        _userdata.TimeoutException,
    ) as e:
      raise ValueError(
          f'Cannot find or access API key in Colab secrets with {secret=}. Make'
          ' sure you have added the API key to Colab secrets and enabled'
          ' access. See'
          ' https://www.alphagenomedocs.com/installation.html#add-api-key-to-secrets'
          ' for more details.'
      ) from e

  raise ValueError(f'Cannot find API key with {secret=}.')
