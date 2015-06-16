#lang racket

(require (only-in racket/class
                  class*))

(require (only-in srfi/26
                  cut))

(provide integrator%)


(define integrator%
  (class object%
         (init universe
               dt
               [force-fields '()])
         (super-new)

         (define _universe universe)
         (define/public (get-universe) _universe)

         (define _dt dt)
         (define/public (get-dt) _dt)

         (define _force-fields (flatten force-fields))
         (define/public (get-force-fields) _force-fields)

         (define (_apply-force-fields paths)
           (let ([total-energy 0.0])
             (for ([ff _force-fields])
               (let-values ([(e gs fcs) (send ff evaluate _universe)])
                 (set! total-energy (+ e total-energy))
                 (for ([idx+g gs])
                   (match-define (cons idx g) idx+g)
                   (send _universe apply-force idx (* -0.5 _dt g)))))
             (send _universe set-potential-energy total-energy)))

         (define/public (step [callback #f])
                        (let ([temperature (send _universe get-temperature)]
                              [paths (send _universe get-paths)])
                          (for-each (cut send <> thermostat temperature _dt)
                                    paths)
                          (_apply-force-fields paths)
                          (for-each (cut send <> propagate _dt)
                                    paths)
                          (_apply-force-fields paths)
                          (for-each (cut send <> thermostat temperature _dt)
                                    paths))
                        (when callback (callback _universe)))))
