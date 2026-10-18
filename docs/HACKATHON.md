# Escape Labs: hackathon submission

**Submitted successfully on October 17, 2026.** Cerebral Valley confirmed the entry and displayed Escape Labs under Your Submission, with Khalifa Ibrahim (Ishaq on the event profile) as the solo participant.

- Track: **Done Means Done**.
- Sponsor category: **Best Use of Nemotron**.
- [Submitted project](https://cerebralvalley.ai/e/crusoe-recursive-agents-hackathon/hackathon/projects).
- [Interactive demonstration](https://escape-lab-three.vercel.app).
- [Submitted Loom video](https://www.loom.com/share/7f7cd52a924a4626af662e5e756726e3).
- [Public implementation and experiment record](https://github.com/Gentle-mann/escape-lab).

Our ambition is an autonomous research lab that anticipates cancer's escape routes, tests its hypotheses, and improves through measured experiments. The submitted six-hour build establishes the research workflow and the matched activation-ranking improvement. Later development is identified separately in the results and commit history.

## Event requirements

Checked against the [official participant guide](https://cerebralvalley.ai/e/crusoe-recursive-agents-hackathon/details) on October 17, 2026.

- Submissions are due **5:00 PM PT, Saturday October 17**.
- Repository must be public; team maximum is four.
- Submit an accessible **one-minute demo video** and add every team member.
- Prepare an approximately **three-minute live demo**, followed by 1–2 minutes of questions.
- Clearly identify code and functionality created during the event. Do not portray pretrained models, published data, or previous work as original contributions. Git publication time alone does not establish when work was built.

## Recommended track and demo

Lead with **Done Means Done**: an agent proposes a bounded experiment, an independent evaluator measures it, rejects unsupported claims, and records the outcome. Crusoe inference, the numerical tool and the verification process should be visible in the demo. A progress/result interface supports this workflow; a dashboard must not be the main feature.

Show the two measured advances with their respective tests: the six-hour numerical ranking improvement on fixed FGFR2 development folds, and the post-submission Nemotron policy loop that improved diagnostic success from 12/48 to 48/48 on generated numerical cases. Nemotron plans investigations and proposes reusable rules; numerical tools produce predictions; the independent evaluator measures each result. An independent biological assay is the next milestone.

Opening line: “An autonomous cancer-research agent that has to earn every claim with a measurement.” Immediately label the project **preclinical research; no patients or treatment recommendations**.

Show one real trace: research evidence, proposed configuration, numerical evaluation, and an accept/reject decision. The AP before/after explains the measured gain: the same classifier had many rounded scores tied at one, and using its raw decision margin recovered ranking information. The complete experiment ledger includes the structure and alternative-planner comparisons.

## Delivery record

The interactive replay, architecture presentation, Loom video, public repository, and submission are complete. Public repository access and signed-out video playback were verified before submission. A three-minute pitch and judge Q&A are prepared. The submitted video displayed approximately 1:12 during review; the builder chose that link with the guide's one-minute request noted.

The demonstration remains a replay of the submitted scientific experiments. Its optional live API check uses synthetic fixtures. The frontend and presentation are delivered separately from this repository; their source remains local under the existing publication arrangement.

The guide prohibits medical-advice generation, basic RAG, Streamlit apps, and dashboard-led projects. Our research framing avoids medical-advice functionality, but this is an interpretation of the guide, not organizer preapproval. Do not add patient-input forms or treatment recommendations to the frontend.

Nemotron and You.com integrations are implemented. NemoClaw/OpenShell were researched but are not deployed in this prototype; do not claim their sponsor integration.
