#include <cmath>
#include <iostream>
#include "Model.h"
#include "MoveRoot.h"
#include "Msg.h"
#include "Node.h"
#include "ParameterTree.h"
#include "RandomVariable.h"
#include "TransitionProbabilities.h"

#undef DEBUG_ROOT_MOVE


MoveRoot::MoveRoot(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]) : Move(r, m, nm) {

    myParameter[0] = t[0];
    myParameter[1] = t[1];
    tuning = log(4.0);
}

void MoveRoot::accept(void) {
    
    numTries++;
    numAccepted++;

    (*myParameter[0]) = (*myParameter[1]);
}

void MoveRoot::reject(void) {

    numTries++;

    (*myParameter[1]) = (*myParameter[0]);
    myParameter[1]->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();
}

void MoveRoot::restore(void) {

}

void MoveRoot::tune(void) {

}

double MoveRoot::update(void) {

    ParameterTree* t = myParameter[1];
    Node* r = t->getRoot();
    std::vector<Node*> rootDescendants = r->getDescendants();
    if (rootDescendants.size() != 2)
        Msg::error("Expecting two descendants of the root node");
    Node* ingroupRoot = NULL;
    Node* outgroupRoot = NULL;
    if (rootDescendants[0]->getIsOutgroup() == true && rootDescendants[1]->getIsOutgroup() == false)
        {
        ingroupRoot = rootDescendants[1];
        outgroupRoot = rootDescendants[0];
        }
    else if (rootDescendants[0]->getIsOutgroup() == false && rootDescendants[1]->getIsOutgroup() == true)
        {
        ingroupRoot = rootDescendants[0];
        outgroupRoot = rootDescendants[1];
        }
    if (ingroupRoot == NULL)
        Msg::error("Could not find ingroup root");
    if (outgroupRoot == NULL)
        Msg::error("Could not find outgroup root");
    std::vector<Node*> ingroupRootDescendants = ingroupRoot->getDescendants();
    std::vector<Node*> outgroupRootDescendants = outgroupRoot->getDescendants();
    if (ingroupRootDescendants.size() != 2)
        Msg::error("Expecting two descendants of the ingroup root node");
    if (outgroupRootDescendants.size() != 0 && outgroupRootDescendants.size() != 2)
        Msg::error("Expecting zero or two descendants of the outgroup root node");


    std::set<Branch*> disallowedBranches;
    disallowedBranches.insert( t->findBranch(r, rootDescendants[0]) );
    disallowedBranches.insert( t->findBranch(r, rootDescendants[1]) );
    disallowedBranches.insert( t->findBranch(ingroupRoot, ingroupRootDescendants[0]) );
    disallowedBranches.insert( t->findBranch(ingroupRoot, ingroupRootDescendants[1]) );
    if (outgroupRootDescendants.size() == 2)
        {
        disallowedBranches.insert( t->findBranch(outgroupRoot, outgroupRootDescendants[0]) );
        disallowedBranches.insert( t->findBranch(outgroupRoot, outgroupRootDescendants[1]) );
        }

    // pick an allowable branch
    bool foundGoodBranch = false;
    Branch* b = NULL;
    do
        {
        b = t->randomBranch();
        std::set<Branch*>::iterator it = disallowedBranches.find(b);
        if (it == disallowedBranches.end())
            foundGoodBranch = true;

        } while(foundGoodBranch == false);
    
    // remove the old root
    Node* d1 = ingroupRootDescendants[0];
    Node* d2 = ingroupRootDescendants[1];
    Branch* b1 = t->findBranch(ingroupRoot, d1);
    Branch* b2 = t->findBranch(ingroupRoot, d2);
    double len = b1->getLength() + b2->getLength();
    t->removeBranch(b2);
    b1->setLength(len);
    b1->setEnds(d1, d2);
    d1->removeNeighbor(ingroupRoot);
    d2->removeNeighbor(ingroupRoot);
    ingroupRoot->removeNeighbor(d1);
    ingroupRoot->removeNeighbor(d2);
    d1->addNeighbor(d2);
    d2->addNeighbor(d1);
    
    // add the new root
    d1 = b->getEnd1();
    d2 = b->getEnd2();
    d1->removeNeighbor(d2);
    d2->removeNeighbor(d1);
    ingroupRoot->addNeighbor(d1);
    ingroupRoot->addNeighbor(d2);
    d1->addNeighbor(ingroupRoot);
    d2->addNeighbor(ingroupRoot);
    b->setEnds(d1, ingroupRoot);
    b2 = t->addBranch(d2, ingroupRoot);
    len = b->getLength();
    double len2 = len * rv->uniformRv();
    b->setLength(len - len2);
    b2->setLength(len2);
    
    // get down pass sequence
    t->initializeDownPassSequence();
    
    // update the flags for conditional likelihoods
    t->flipAllActiveConditionalLikelihoods();
    t->updateAllConditionalLikelihoods(true);

    // update the transition probability for the branch
    t->updateAllTransitionProbabilities(true);
    model->updateTransitionProbabilities();

    return 0.0;
}

std::vector<double> MoveRoot::values(void) {

    std::vector<double> vals;
    vals.push_back( 0.0 );
    return vals;
}
