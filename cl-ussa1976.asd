;;;; -*- Mode: Lisp -*-

(in-package :cl-user)

(defpackage cl-ussa1976-system
  (:use :common-lisp :asdf))

(in-package :cl-ussa1976-system)

(defsystem cl-ussa1976
  :description "US Standard Atmosphere, 1976"
  :version "0.1"
  :author "Mayes Mullins <mmullins@mullinsenterprises.ca>"
  :licence "MIT"
  :depends-on ()
  :components ((:file "packages")
               (:file "cl-ussa1976" :depends-on ("packages"))))
