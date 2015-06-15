#lang racket

(require
  (only-in srfi/26
           cut)
  (only-in srfi/43
           [vector-map! vector-map-enum!]
           vector-reverse-copy!)
  (only-in (planet williams/science/fft)
           fft-complex-forward
           fft-complex-inverse))

(provide dct
         idct)


(define (vector-update! v n f)
  (vector-set! v n
               (f (vector-ref v n))))

(define (map-etotheipi! v x)
  (vector-map-enum!
    (lambda (n y)
      (* (exp (/ (* pi x n)
                 (vector-length v)))
         y))
    v))

;; Type-II discrete cosine transform (DCT-II) in 1 dimension. Also known as
;; "the DCT".
(define (dct xs)
  (define N (vector-length xs))
  (define xs2 (make-vector (* 2 N)))
  (vector-copy! xs2 0 xs)
  (vector-reverse-copy! xs2 N xs)
  (fft-complex-forward xs2)
  (define Xs (vector-take xs2 N))
  (map-etotheipi! Xs -0.5i)
  (vector-map! real-part Xs)
  (vector-update! Xs 0 (cut / <> (sqrt 2)))
  (vector-map! (cut / <> (sqrt (* 2 N))) Xs))

;; Type-III discrete cosine transform (DCT-III) in 1 dimension. Also known as
;; "the inverse DCT".
(define (idct Xs)
  (define N (vector-length Xs))
  (define Xs2 (vector-copy Xs))
  (map-etotheipi! Xs2 +0.5i)
  (vector-update! Xs2 0 (cut * <> (sqrt 2)))
  (define Xs3 (make-vector (* 2 N)))
  (vector-copy! Xs3 0 Xs2)
  (vector-set! Xs3 N 0)
  (vector-reverse-copy! Xs3 (add1 N)
                        (vector-map conjugate Xs2) 1)
  (vector-map! (cut * <> (sqrt (* 2 N))) Xs3)
  (fft-complex-inverse Xs3)
  (vector-map real-part (vector-take Xs3 N)))
