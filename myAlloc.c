/**
* Writer: Xiaolin.liu
* xiaolin.liu@mail.bnu.edu.cn
*
* This module contains basic functions for  calculation.
* Functions list:
* Kernel: 
* 20xx.xx.xx, LOC
**/

#include "myAlloc.h"
#include <pthread.h>
static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

typedef struct tagAllocNodeList
{
    void *addr;
    size_t size;
    const char *file;
    int line;
    struct tagAllocNodeList *next;
} AllocNodeList;

// A node list to record all memory allocated
static AllocNodeList *alloc_node_head = NULL;

// Check allocation linked node list and free all nodes remain.
void CheckMemoryLeak()
{
    int leak;
    AllocNodeList *i = NULL, *this = NULL;
    i = alloc_node_head;
    while (i)
    {
        this = i;
        if (this->addr != NULL)
        {
            fprintf(stderr, "%s_ %p: %zu bytes (%s:%d)\n", __func__, this->addr, this->size, this->file, this->line);
            free(this->addr);
        }
        free(this);
        this = NULL;
        i = i->next;
    }
    return;
}

// When a new memory block was allocated, record it in the list.
static void AddNodeToList(AllocNodeList *new)
{
    if (!alloc_node_head)
    {
        alloc_node_head = new;
        return;
    }
    AllocNodeList *i = NULL;
    i=alloc_node_head;
    while(i->next)
    {
        i = i->next;
    }
    i->next = new;
    return;
}

// When a recorded memory block was freed, remove it from the list.
static int RemoveNodeFromList(void *addr)
{
    if (!alloc_node_head)
    {
        fprintf(stderr, "%p is not in the node list\n", addr);
        return 0;
    }
    AllocNodeList *now = NULL, *prev = NULL;
    INT icount = 0;
    now=alloc_node_head;
    while(now != NULL)
    {
        icount++;
        if (now->addr == addr)
        {
            break;
        }
        prev = now;
        now = now->next;
    }

    if (now == NULL)
    {
        fprintf(stderr, "%p is not in the node list\n", addr);
        return 0;
    }
    else if (icount == 1)
    {
        alloc_node_head = now->next;
        free(now);
        now = NULL;
        return 1;
    }

    prev->next = now->next;
    free(now);
    now = NULL;
    return 1;
}

static void *PushAlloc(void *p, size_t n, const char *file, int line)
{
    AllocNodeList *newNode;
    if (!(newNode = malloc(sizeof(*newNode))))
    {
        return NULL;
    }
    pthread_mutex_lock(&mut);
    newNode->addr = p;
    newNode->size = n;
    newNode->file = file;
    newNode->line = line;
    newNode->next = NULL;
    AddNodeToList(newNode);
    pthread_mutex_unlock(&mut);
    return p;
}

void *myMallocLong(size_t n, const char *file, int line)
{
    void *p = NULL;
    void *q = NULL;
    p = malloc(n);
    q = PushAlloc(p, n, file, line);
    if (!q)
    {
        fprintf(stderr, "[%s] [%d]Failed to alloc memory\n", file, line);
        if (p)
            free(p);
    }
    return q;
}

void *myCallocLong(size_t m, size_t n, const char *file, int line)
{
    size_t sz;
    void *p;
    void *q;
    sz = m * n;
    p = malloc(sz);
    q = PushAlloc(p, sz, file, line);
    if (!q)
    {
        fprintf(stderr, "[%s] [%d]Failed to alloc memory\n", file, line);
        if (p)
            free(p);
    }
    return q ? memset(q, 0, sz) : NULL;
}

void myFree(void *p)
{
    if (!p)
        return;
    RemoveNodeFromList(p);
    free(p);
    return;
}