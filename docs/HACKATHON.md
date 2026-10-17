# Hackathon delivery checklist

Checked against the [official participant guide](https://cerebralvalley.ai/e/crusoe-recursive-agents-hackathon/details) on October 17, 2026.

- Submissions are due **5:00 PM PT, Saturday October 17**.
- Repository must be public; team maximum is four.
- Submit an accessible **one-minute demo video** and add every team member.
- Prepare an approximately **three-minute live demo**, followed by 1–2 minutes of questions.
- Clearly identify code and functionality created during the event. Do not portray pretrained models, published data, or previous work as original contributions. Git publication time alone does not establish when work was built.

## Recommended track and demo

Lead with **Done Means Done**: an agent proposes a bounded experiment, an independent evaluator measures it, rejects unsupported claims, and records the outcome. Crusoe inference, the numerical tool and the verification process should be visible in the demo. A progress/result interface supports this workflow; a dashboard must not be the main feature.

Treat Version N+1 as an implemented mechanism, not proven unseen-task improvement. The current AP repair was developed during analysis, not discovered by Nemotron or GLM. Reused development gains cannot satisfy the track's unseen-task evidence by themselves.

Opening line: “An autonomous cancer-research agent that has to earn every claim with a measurement.” Immediately label the project **preclinical research; no patients or treatment recommendations**.

Show one real trace: research evidence, proposed configuration, numerical evaluation, accept/reject decision, and recorded limitations. The AP repair can be a separate before/after inspection: the same classifier had many rounded scores tied at one, and using its raw decision margin recovered ranking information. Show the unsuccessful structure and LLM-swap experiments in the ledger rather than hiding them.

## Remaining deliverables

1. Build a presentation interface around the actual agent/evaluator workflow, with provenance and clear live versus recorded state.
2. Rehearse one reliable run using development data, with a recorded fallback if inference is slow. No new scientific parameter sweep is needed for the demo.
3. Prepare the one-minute video and a three-minute talk track; verify video access from a signed-out session.
4. Add a concrete contribution statement listing what was built at the event and what is reused.
5. Recheck public repository access, correct experiment claims, team members, and links; submit at the [official form](https://cerebralvalley.ai/e/crusoe-recursive-agents-hackathon/hackathon/submit).

The guide prohibits medical-advice generation, basic RAG, Streamlit apps, and dashboard-led projects. Our research framing avoids medical-advice functionality, but this is an interpretation of the guide, not organizer preapproval. Do not add patient-input forms or treatment recommendations to the frontend.

Nemotron and You.com integrations are implemented. NemoClaw/OpenShell were researched but are not deployed in this prototype; do not claim their sponsor integration.
