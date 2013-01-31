#ifndef QUEUE_H_
#define QUEUE_H_

#include <exception> // general exception
#include <stdexcept> // out_of_range
#include <vector>

/* exceptions */
class QueueEmptyException : public std::out_of_range {
    public:
        QueueEmptyException(const std::string &message) : std::out_of_range(message) {;}
};

/**
 * It's a queue.
 */

template <class Type>
struct node {
    node* next;
    Type data;
};

template <class Type>
class Queue {
    private:
        node<Type> *front, *back;
        int currentSize;
    public:
        Queue() {
            front = back = NULL;
            currentSize = 0;
        }

        void push(Type item) {
            /*  Empty queue */
            if (front == NULL) {
                front = new node<Type>;
                back = front;
            } else {
                back->next = new node<Type>;
                back = back->next;
                back->next = NULL;
            }
            back->data = item;
            currentSize++;
        }

        Type pop() {
            if (!currentSize) {
                throw QueueEmptyException("There isn't anything left in the queue ;_;");
            }
            Type ret = front->data;
            node<Type>* temp = front;
            front = front->next;
            delete temp;
            currentSize--;
            return ret;
        }

        Type peek() {
            if (front == NULL) {
                throw QueueEmptyException("There isn't anything left in the queue ;_;");
            }
            return front->data;
        }

        int size() {
            return currentSize;
        }

};

#ifdef MATLAB_MEX_FILE
#include "mex.h"

void retrieve_queue(const mxArray* matptr, Queue<std::vector<double> >* &queue){
    // retrieve pointer from the MX form
    double* pointer0 = mxGetPr(matptr);
    if (!pointer0)
        mexErrMsgTxt("vararg{1} must be a valid queue pointer\n");

    // convert it to "long" datatype (make this void*)
    long pointer1 = (long) pointer0[0];
    // convert it to "queue"
    queue = (Queue<std::vector<double> >*) pointer1;
    if (!queue)
        mexErrMsgTxt("vararg{1} must be a valid queue pointer\n");
}
#endif

#endif /*QUEUE_H_*/
