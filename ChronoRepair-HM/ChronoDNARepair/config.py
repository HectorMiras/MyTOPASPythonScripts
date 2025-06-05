#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np

class Config:
    def __init__(self, filename):
        config_dict = self.read_config_file(filename)
        self.__dict__ = {**self.__dict__, **config_dict}

        if self.seed == -1:
            self.seed = np.random.randint(0, 1000000)

        self.working_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        print('Working directory: ' + self.working_directory)

    def read_config_file(self, filename):
        config_dict = {}
        filename = filename + '.txt'

        with open(filename, 'r') as f:
            for line in f:
                # Remove comments and strip leading/trailing whitespace
                line = line.split('#')[0].strip()
                # Skip blank lines
                if line == '':
                    continue
                # Split the line into key/value pairs
                key, value = line.split(': ')
                # Check if the value is a number (integer or float)
                if value.replace('.', '', 1).isdigit() or value.lstrip('-').replace('.', '', 1).isdigit():
                    if '.' in value:
                        config_dict[key] = float(value)
                    else:
                        config_dict[key] = int(value)
                # Check if the value is a boolean
                elif value.lower() == 'true' or value.lower() == 'false':
                    config_dict[key] = bool(value.lower() == 'true')
                # Otherwise, treat the value as a string
                else:
                    config_dict[key] = value
        return config_dict