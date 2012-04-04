;;;; -*- Mode: Lisp -*-

(in-package :cl-user)

(defpackage test-cl-ussa1976-system
  (:use :common-lisp :asdf))

(in-package :test-cl-ussa1976-system)

(defsystem test-cl-ussa1976
    :description "Test support for cl-ussa1976.  Functions to generate tables like those in the standard."
    :version "0.1"
    :author "Mayes Mullins <mmullins@mullinsenterprises.ca>"
    :licence "MIT"
    :depends-on (:cl-ussa1976)
    :components ((:file "packages")
                 (:file "test-cl-ussa1976" :depends-on ("packages"))))
