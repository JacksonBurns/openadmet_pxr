todo:
 - make a chemical space projection that follows the below trajectory
 - check distribution of pEC50 for PXR for compounds inside the counter screen and those not to validate if step 2 is a good idea (and do the above)
 - we have uncertainty bounds - incorporate into loss function?

current plan:

ensemble of the below models:
A. Chemprop
1. pre-train chemprop on the primary screen with concentration as an input descriptor target is log2_fc_estimate
2. fine tune (2) on compounds put through pxr and counter (multitask)

B. CheMeleon
Direct training on the pec50 task (multitask as well)

C. Physciochemical Random Forest

From the CheMeleon paper, on pEC50 only since multitask doesn't help here

later considerations - can add uncertainty to the loss function with the known uncertainty values, possibly

note for pretrained model:
Avoid data leaks by dropping members from earlier datasets that appear in later datasets.
