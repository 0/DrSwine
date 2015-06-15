#lang racket

(require (only-in racket/class
                  class*))

(require (only-in srfi/26
                  cut))

(provide integrator%)


(define integrator%
  (class object%
         (init universe
               dt)
         (super-new)

         (define _universe universe)
         (define/public (get-universe) _universe)

         (define _dt dt)
         (define/public (get-dt) _dt)

         (define/public (step [callback #f])
                        (let ([temperature (send _universe get-temperature)]
                              [paths (send _universe get-paths)])
                          (for-each (cut send <> thermostat temperature _dt)
                                    paths)
                          (for-each (cut send <> propagate _dt)
                                    paths)
                          (for-each (cut send <> thermostat temperature _dt)
                                    paths))
                        (when callback (callback _universe)))))
