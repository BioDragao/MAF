#lang racket/base
(require racket/date txexpr pollen/setup)
(provide (all-defined-out))
 
(module setup racket/base
  (provide (all-defined-out))
  (define poly-targets '(html )))
 
(define (get-date)
  (date->string (current-date)))
 