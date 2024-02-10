#!/bin/bash
timeout --preserve-status 60s fpm run Randomized
exit $?
