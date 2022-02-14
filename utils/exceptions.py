"""
This file holds all ammonia-plant exceptions, organised by exception id
"""


class AmmoniaPlantException(Exception):
    """
    The base class for all ammonia-plant exceptions

    To define a new exception, inherrit from this and override the
    _class_message
    """
    class_message = 'An unknown exception occurred.'
    error_code = 1

    def __init__(self, message=''):
        super().__init__(message)
        self._obj_message = message
        self._error_string = None

    def __str__(self):
        if self._obj_message != '':
            self._error_string = '{}\nDetails: {}'.format(self.class_message,
                                                          self._obj_message)
        else:
            self._error_string = self.class_message

        return self._error_string.strip()

class OptionsError(AmmoniaPlantException):
    """
    Indicates an error during processing options.
    """
    class_message = 'Failed to process options.'
    error_code = 2