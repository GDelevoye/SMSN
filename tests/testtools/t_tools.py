from itertools import zip_longest

### Helper functions for the tests that are not really pytest fixtures

def zip_equal(*iterables):
    """Enforces that the iterables of a zip must have the same length
    Taken from SO https://stackoverflow.com/questions/32954486/zip-iterators-asserting-for-equal-length-in-python

    >>> liste1 = [0,1,2,3]
    >>> liste2 = [4,5,6,7]
    >>> liste3 = [8,9,10]
    >>> # https://stackoverflow.com/questions/12592/can-you-check-that-an-exception-is-thrown-with-doctest-in-python
    >>> [(x,y) for (x,y) in zip_equal(liste1,liste2)]
    [(0, 4), (1, 5), (2, 6), (3, 7)]
    >>> [(x,y) for (x,y) in zip_equal(liste1,liste3)]
    Traceback (most recent call last):
     ...
    ValueError: Iterables have different lengths
    """
    sentinel = object()
    for combo in zip_longest(*iterables, fillvalue=sentinel):
        if sentinel in combo:
            raise ValueError('Iterables have different lengths')
        yield combo

###
