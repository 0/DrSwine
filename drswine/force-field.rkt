#lang racket

(require (only-in racket/class
                  class*))

(require (only-in srfi/26
                  cut))

(provide apply-ff-to-indices
         apply-ff-to-paths
         ff-harmonic-distance
         ff-harmonic-trap
         force-field%)


(define force-field%
  (class* object% (printable<%>)
          (init f
                affected-indices
                [scaling 1.0]
                [name #f])
          (super-new)

          (define _f f)

          (define _affected-indices affected-indices)
          (define/public (get-affected-indices) _affected-indices)

          (define _scaling scaling)
          (define/public (get-scaling) _scaling)
          (define/public (set-scaling new-scaling)
                         (set! _scaling new-scaling))

          (define _name name)
          (define/public (get-name) _name)

          (define/public (evaluate universe)
                         (let-values ([(e gs fcs) (_f universe)])
                           (values (* _scaling e)
                                   (for/list ([g gs])
                                             (cons (car g) (* _scaling (cdr g))))
                                   (for/list ([fc fcs])
                                             (cons (car fc) (* _scaling (cdr fc)))))))

          (define/public (repr) (format "#<~a ~a ~a>"
                                        _name
                                        (set->list _affected-indices)
                                        _scaling))

          (define/public (custom-print p depth) (display (repr) p))
          (define/public (custom-write p)       (write (repr) p))
          (define/public (custom-display p)     (display (repr) p))))


(define (apply-ff-to-indices ff indices
                             #:scale [scaling 1.0]
                             #:name [name #f])
  (new force-field%
       [f (apply ff indices)]
       [affected-indices (apply set indices)]
       [scaling scaling]
       [name name]))

(define (apply-ff-to-paths ff path
                           #:scale-all [scale-all #f]
                           #:scale-ends [scale-ends #f]
                           #:only-ends [only-ends #f]
                           #:name [name #f]
                           . more-paths)
  (let* ([paths (cons path more-paths)]
         [lengths (for/list ([p paths])
                            (send p get-num-beads))]
         [indices (for/vector ([p paths])
                              (send p get-bead-indices))]
         [len (car lengths)]
         [scale-all (if (eq? 'links scale-all) (/ 1.0 len) scale-all)]
         [scale-ends (if (eq? 'pigs scale-ends) 0.5 scale-ends)]
         [scaling (make-vector len
                               (if scale-all scale-all 1.0))]
         [use-indices (if only-ends
                        (vector 0 (sub1 len))
                        (build-vector len identity))])
    (unless (or (null? (cdr lengths))
                (apply = lengths))
      (error "paths must be of equal lengths"))
    (when scale-ends
      (vector-set! scaling 0
                   (* scale-ends
                      (vector-ref scaling 0)))
      (vector-set! scaling (sub1 len)
                   (* scale-ends
                      (vector-ref scaling (sub1 len)))))
    (for/list ([i use-indices])
              (let ([the-indices
                      (vector->list
                        (vector-map (cut vector-ref <> i) indices))])
                (new force-field%
                     [f (apply ff the-indices)]
                     [affected-indices (apply set the-indices)]
                     [scaling (vector-ref scaling i)]
                     [name name])))))


; Harmonic central trap potential for a single bead.
(define (((ff-harmonic-trap k pos) idx) u)
  (let* ([q (send u get-position idx)]
         [dq (- q pos)]
         [energy (* 0.5 k (sqr dq))]
         [gradients (list (cons idx (* k dq)))]
         [fcs (list (cons idx k))])
    (values energy gradients fcs)))


; Harmonic distance potential for a pair of beads.
(define (((ff-harmonic-distance k) idx1 idx2) u)
  (let* ([q1 (send u get-position idx1)]
         [q2 (send u get-position idx2)]
         [dq (- q1 q2)]
         [energy (* 0.5 k (sqr dq))]
         [gradients (list (cons idx1 (* +1.0 k dq))
                          (cons idx2 (* -1.0 k dq)))]
         [fcs (list (cons idx1 k)
                    (cons idx2 k))])
    (values energy gradients fcs)))
